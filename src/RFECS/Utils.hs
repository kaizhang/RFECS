{-# LANGUAGE OverloadedStrings #-}

module RFECS.Utils
    ( withTmpDir
    , withTmpFile
    , readLoops
    , readChrSize
    ) where

import qualified Data.ByteString.Char8 as B
import Bio.Data.Bed
import Bio.Utils.Misc
import           Control.Exception (bracket)
import qualified Data.Text         as T
import           Shelly            hiding (FilePath, withTmpDir)

withTmpFile :: FilePath -> (FilePath -> IO a) -> IO a
withTmpFile dir = bracket create delete
  where
    create = shelly $ fmap (T.unpack . head . T.lines) $ silently $
        run "mktemp" [T.pack $ dir ++ "/RFECS_tmp_XXXXXXXX"]
    delete = shelly . rm . fromText . T.pack

withTmpDir :: FilePath -> (FilePath -> IO a) -> IO a
withTmpDir dir = bracket create delete
  where
    create = shelly $ fmap (T.unpack . head . T.lines) $ silently $
        run "mktemp" ["-d", T.pack $ dir ++ "/RFECS_tmp_XXXXXXXX"]
    delete = shelly . rm_rf . fromText . T.pack

readLoops :: FilePath -> IO [(BED3, BED3)]
readLoops fl = (map toLoop . tail . B.lines) <$> B.readFile fl
{-# INLINE readLoops #-}

toLoop :: B.ByteString -> (BED3, BED3)
toLoop x =
    let xs = B.split '\t' x
    in ( asBed (xs !! 0) (readInt $ xs !! 1) (readInt $ xs !! 2)
       , asBed (xs !! 3) (readInt $ xs !! 4) (readInt $ xs !! 5)
       )
{-# INLINE toLoop #-}

readChrSize :: FilePath -> IO [(B.ByteString, Int)]
readChrSize fl = do
    c <- B.readFile fl
    return $ map f $ B.lines c
  where
    f x = let [chr,size] = B.words x
          in (chr, readInt size)
{-# INLINE readChrSize #-}
