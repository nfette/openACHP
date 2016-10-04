FOR %%i IN (00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19) DO montage convergeme2.*.%%i.png -geometry 400x300 montage2.%%i.png
convert montage2.* anewfile2.pdf
