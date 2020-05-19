{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
module Main where

import           Bio.ChIPSeq                       (rpkmBinBed)
import           Bio.Data.Bed
import           Bio.Utils.Functions               (slideAverage)
import           Bio.Utils.Misc                    (readDouble, readInt)
import           Conduit
import           Control.Arrow                     (second)
import           Control.Lens                      ((.~), (^.))
import           Control.Monad
import           Control.Monad.Base                (liftBase)
import           Control.Monad.Morph               (hoist)
import qualified Data.ByteString.Char8             as B
import           Data.Conduit.Zlib                 (ungzip)
import           Data.Default
import           Data.Double.Conversion.ByteString (toFixed)
import           Data.List
import qualified Data.Matrix.Unboxed               as MU
import           Data.Maybe
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import           Data.Version                      (showVersion)
import           Data.Yaml                         (decodeFileEither,
                                                    prettyPrintParseException)
import           Options.Applicative
import           Paths_RFECS                       (version)
import           Shelly                            hiding (FilePath, withTmpDir)
import           System.IO
import           Text.Printf                       (printf)

import           RFECS.Constants
import           RFECS.Type                        hiding (info)
import           RFECS.Utils

data Options = Options
    { input   :: FilePath
    , output  :: FilePath
    , chrSize :: String
    , cmd     :: FilePath
    , chrsets :: String
    , tmp :: FilePath
    } deriving (Show, Read)

parser :: Parser Options
parser = Options
     <$> argument str (metavar "INPUT")
     <*> argument str (metavar "OUTPUT")
     <*> strOption
           ( long "genome"
          <> short 'g'
          <> help "A file that lists all chromosomes' sizes. Built-in genome: hg19 and mm10."
          <> metavar "GENOME" )
     <*> strOption
           ( long "path"
          <> short 'p'
          <> help "The path to the core RFECS matlab scripts which are provided separately."
          <> metavar "RFECS_PATH" )
     <*> strOption
           ( long "chrom"
          <> help "Specify which chromosomes to be analyzed, example: 'chr1,chr2'. Use 'all' for whole genome. Default: all."
          <> value "all"
          <> metavar "CHROMOSOME" )
     <*> strOption
           ( long "tmp"
          <> help "Specify the directory for storing temparary files. Default: current directory."
          <> value "./"
          <> metavar "TMP_DIR" )

defaultMain :: Options -> IO ()
defaultMain (Options inFl output chrsize rfecs chrsets tmp) = do
    genome <- case chrsize of
        "hg19" -> return hg19ChrSize
        "mm10" -> return mm10ChrSize
        _ -> readChrSize chrsize

    let chr | chrsets == "all" = genome
            | otherwise = map (\x -> (x, fromJust $ lookup x genome)) $
                B.split ',' $ B.pack chrsets

    r <- decodeFileEither inFl
    case r of
        Left e -> error $ prettyPrintParseException e
        Right inputs -> withTmpDir tmp $ \tmpDir -> do
            let chrBed = map (\(chr,e) -> BED3 chr 0 e) chr
            rcs <- readCount chrBed inputs
            writeReadCount tmpDir $ zip (fst $ unzip chr) rcs
            rfecs_matlab "matlab" rfecs (rfecs ++ "/model_human_3marks.mat")
                tmpDir tmpDir
            fls <- shelly $ lsT $ fromText $ T.pack tmpDir
            enhancers <- forM (filter ("_enhancers.txt" `T.isSuffixOf`) fls) $ \fl -> do
                c <- B.readFile $ T.unpack fl
                return $ map (f . B.split '\t') $ B.lines c
            writeBed' output $ concat enhancers
  where
    f [a,b,c] = BED a (readInt b - 1000) (readInt b + 1000) Nothing
        (Just $ readDouble c) Nothing


writeReadCount :: FilePath -> [(B.ByteString, MU.Matrix Float)] -> IO ()
writeReadCount dir xs = forM_ xs $ \(chr, rc) -> do
    shelly $ mkdir_p $ fromText $ T.pack dir
    withFile (dir ++ "/" ++ B.unpack chr) WriteMode $ \h -> do
        let n = MU.rows rc
        forM_ [0..n-1] $ \i -> B.hPutStrLn h $ B.intercalate "\t" $
            map (toFixed 6 . realToFrac) $ U.toList $ MU.takeRow rc i
{-# INLINE writeReadCount #-}

readCount :: [BED3]  -- ^ regions
          -> [Experiment]
          -> IO [MU.Matrix Float]  -- ^ each matrix in the list stores
                                   -- read count of one chromosome across different experiment
readCount chr es = do
    rpkms <- forM es $ \exp -> do
        let [fl] = exp^.files
        (v, n) <- runResourceT $ readBedFromFile fl $$ hoist liftBase (rpkmBinBed 100 chr)
        let factor = fromIntegral n / 1e7
            v' :: [U.Vector Float]
            v' | exp^.target == "Input" = map (slideAverage 2 . U.map (/ factor)) v
               | otherwise = map (U.map (/ factor)) v
        return (exp^.eid, v')

    let rs = flip map (filter f rpkms) $ \(id', v) ->
            let e = fromJust $ lookup id' expMap
            in case e ^. control of
                Just c -> (e^.target, zipWith (U.zipWith (-)) v $ fromJust $ lookup c rpkms)
                _ -> (e^.target, v)
        h3k27ac = fromJust $ lookup "H3K27ac" rs
        h3k4me1 = fromJust $ lookup "H3K4me1" rs
        h3k4me3 = fromJust $ lookup "H3K4me3" rs
    return $ zipWith3 (\a b c -> MU.fromColumns [a,b,c]) h3k27ac h3k4me1 h3k4me3
  where
    expMap = map (\x -> (x^.eid, x)) es
    f (id',_) = (fromJust (lookup id' expMap) ^. target) /= "Input"
    readBedFromFile fl = case fl^.format of
        Bed -> readBed $ fl^.location
        BedGZip -> sourceFile (fl^.location) =$= ungzip =$=
            linesUnboundedAsciiC =$= mapC fromLine
        _ -> error "Unknown input format"
{-# INLINE readCount #-}

rfecs_matlab :: FilePath -> FilePath -> FilePath -> FilePath -> FilePath -> IO ()
rfecs_matlab matlab_path rfecs_path model inputDir output = do
    chrset <- fmap (intercalate "," . map ((\x -> "'" ++ x ++ "'") . T.unpack .
        snd . T.breakOnEnd "/")) $ shelly $ lsT $ fromText $ T.pack inputDir
    let script = printf ( "addpath('%s')," ++
            "parpool(2)," ++
            "extn ='_enhancers';," ++
            "chr_set={%s};," ++
            "load %s," ++
            "forest_predict_par(forest_p300_tssbg_all,3,65,20,[1:3],extn,20,'%s/','%s/',chr_set);," ++
            "parpool close," ++
            "exit" ) rfecs_path chrset model inputDir output
    shelly $ run_ (fromText $ T.pack matlab_path)
        ["-nodisplay", "-nosplash", "-nodesktop", "-r", T.pack script]
  where
{-# INLINE rfecs_matlab #-}

main :: IO ()
main = execParser opts >>= defaultMain
  where
    opts = info (helper <*> parser)
            ( fullDesc
           <> header (printf "RFECS-v%s" (showVersion version)) )
