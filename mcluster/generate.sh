make clean 
make mcluster
#./mcluster -N 10000 -B 1000 -s 12345 -u 1 -C 5 -o CometsProject



unit=0 # 0 --> Nbody unit;  1 --> astrophysics unit
Ns=128000
Nc=$Ns
#B=51200
B=0
T=1000
NAME=Ns${Ns}_B${B}_comets${Nc}_timeNB_$T


rm *.info *.input *.10


./mcluster -N $Ns -B $B -m 0.08 -m 10 -s 12345 -t 3 -u $unit -C 5 -T $T -o $NAME > log.gen
var=`awk '$1=="scalingInfo" {print $4, $7}' log.gen`
echo $var
#python aei.py $var $NAME


Q=0.9
./mcluster -N $Nc  -m 0.08 -m 10 -Q $Q -s 54321 -t 3 -u $unit -C 5 -o ${NAME}COMs &> /dev/null
awk '{print "0.0000000000000001", $2, $3, $4, $5, $6, $7}' *COMs.dat.10 > dat.10
rm ${NAME}COMs.*


cat ${NAME}.dat.10 >> dat.10


rm -rf ../nbody && mkdir ../nbody
cp dat.10 *.input start.lsf start.sh ../nbody

