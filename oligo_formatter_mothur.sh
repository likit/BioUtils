echo "forward " `grep -v \# ${1} | head -1 | cut -f 3` | tr " " "\t" > forward.txt

grep -v \# ${1} | cut -f 2,4 > tmp1

for n in $(cat tmp1)
do
        echo "barcode" >> barcode.txt;
done

paste barcode.txt tmp1 > tmp2
cat forward.txt tmp2 > ${1}.oligo

rm barcode.txt forward.txt tmp1 tmp2
