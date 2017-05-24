make clean 
make mcluster
#./mcluster -N 10000 -B 1000 -s 12345 -u 1 -C 5 -o CometsProject



unit=0 # 0 --> Nbody unit;  1 --> astrophysics unit
N=128000 # Ns: total number of particles or stars(for ffc)
B=12800
Nnomass=128000
T=100
NAME=Nstar${N}_B${B}_comets${Nnomass}_timeNB_$T


rm *.info *.input *.10


./mcluster -N $N -n $Nnomass -B $B -m 0.08 -m 10 -s 12345 -t 3 -u $unit -C 5 -T $T -o $NAME > log.gen
var=`awk '$1=="scalingInfo" {print $4, $7}' log.gen`
echo $var
#python aei.py $var $NAME

echo 'first step finished!'

#Q=0.5
./mcluster -N $Nnomass  -m 0.08 -m 10  -s 54321 -t 3 -u $unit -C 5 -o ${NAME}COMs &> /dev/null
exit
awk '{print "0.0000000000000001", $2, $3, $4, $5, $6, $7}' *COMs.dat.10 > dat.10
#rm ${NAME}COMs.*


echo 'second step finished!'

cat ${NAME}.dat.10 >> dat.10


rm -rf ../nbody && mkdir ../nbody
cp dat.10 *.input start.lsf start.sh ../nbody

