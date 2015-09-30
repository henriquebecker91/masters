import Data.List (sort{-, sortBy-})
import Data.String (fromString)

merge :: [(Int, Int)] -> [(Int, Int)] -> [(Int, Int)]
merge [] ys = ys
merge xs [] = xs
merge xs@(i@(w, p):xss) ys@(i'@(w', p'):yss) =
	if w < w' then i : (merge xss ys)
	else if w > w' then i' : (merge xs yss)
	else (w, max p p') : (merge xss yss)

addtest :: (Int, Int) -> [(Int, Int)] -> Int -> [(Int, Int)]
addtest i@(w, p) [] c = []
addtest i@(w, p) ((w', p'):xss) c =
	if w+w' > c then []
	else (w+w', p+p'):(addtest i xss c)

filter_d :: [(Int, Int)] -> Int -> [(Int, Int)]
filter_d [] _ = []
filter_d (i@(w, p):xss) max =
	if p > max then i : (filter_d xss p)
	else filter_d xss max

{-weight_non_decreasing (a1, b1) (a2, b2)
	| a1 > a2 = GT
	| a1 < a2 = LT
	| a1 == a2 = compare b1 b2-}

s :: Int -> [(Int, Int)] -> Int -> [(Int, Int)]
s 0 _ _ = []
s k ys c = x
	where x = filter_d (merge (s (k-1) ys c) (addtest (ys !! (k-1)) ((0,0):x) c)) 0

s' str = (show $ last $ s n ys' c) ++ "\n"
	where	ls = lines str 
		(n, c) = (read $ ls !! 0, read $ ls !! 1)
		ys = map (\(w:p:[]) -> (read w, read p)) $ map words $ drop 2 ls
		ys' = sort ys

main = interact s'

