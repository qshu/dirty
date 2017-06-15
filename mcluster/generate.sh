make clean 
make mcluster
#./mcluster -N 10000 -B 1000 -s 12345 -u 1 -C 5 -o CometsProject



unit=0 # 0 --> Nbody unit;  1 --> astrophysics unit
N=128000 # Ns: total number of particles or stars(for ffc)
B=0
Nnomass=128000
T=100
NAME=Nstar${N}_B${B}_comets${Nnomass}_timeNB_$T


rm *.info *.input *.10


./mcluster -N $N -n $Nnomass -B $B -m 0.08 -m 10 -s 12345 -t 3 -u $unit -C 5 -T $T -o $NAME > log.gen
var=`awk '$1=="scalingInfo" {print $4, $7}' log.gen`
echo $var
#python aei.py $var $NAME

echo 'first step finished!'

./mcluster -N $Nnomass -n $Nnomass -m 0.08 -m 10  -s 54321 -t 3 -u $unit -C 5 -T $T  > test.log.gen
#./mcluster -N 128000 -n 128000 -B 12800 -m 0.08 -m 10 -s 54321 -t 3 -u 0 -C 5 -T 100 > test.log.gen
#exit
awk 'BEGIN {k = sqrt(0.5/0.5)} { print $1 * 0.000000001, $2, $3, $4, $5 * k, $6 * k, $7 * k}' test.dat.10 > dat.10
#awk 'BEGIN {k = sqrt(0.5/0.5)} { printf("%20.18f %20.18f %20.18f %20.18f %20.18f %20.18f %20.18f\n" $1 * 0.000000001, $2, $3, $4, $5 * k, $6 * k, $7 * k)}' test.dat.10 > dat.10
#rm ${NAME}COMs.*
rm test.*


echo 'second step finished!'

cat ${NAME}.dat.10 >> dat.10


rm -rf ../nbody && mkdir ../nbody
mv dat.10 *.input ../nbody/
#cp dat.10 *.input start.lsf start.sh ../nbody/

