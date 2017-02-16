make clean 
make mcluster
#./mcluster -N 10000 -B 1000 -s 12345 -u 1 -C 5 -o CometsProject



unit=0 # 0 --> Nbody unit;  1 --> astrophysics unit
N=12800
#B=51200
B=0
NAME=N12.8k_b0_comets0k

rm *.info *.input *.10


./mcluster -N $N -B $B -m 0.08 -m 10 -s 12345 -t 3 -u $unit -C 5 -o $NAME > log.gen
var=`awk '$1=="scalingInfo" {print $4, $7}' log.gen`
echo $var
#python aei.py $var $NAME



#./mcluster -N $N -B $B -m 0.08 -m 10 -s 54321 -u $unit -C 5 -o ${NAME}COMs &> /dev/null
#awk '{print "0.0000000000000001", $2, $3, $4, $5, $6, $7}' *COMs.dat.10 > dat.10
#rm ${NAME}COMs.*


cat ${NAME}.dat.10 >> dat.10


rm -rf ../nbody && mkdir ../nbody
cp dat.10 *.input ../nbody/
