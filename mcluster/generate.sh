make clean 
make mcluster
#./mcluster -N 10000 -B 1000 -s 12345 -u 1 -C 5 -o CometsProject



unit=0 # 0 --> Nbody unit;  1 --> astrophysics unit
Nstar=12800 # Ns: total number of particles or stars(for ffc)
B=0
Nnomass=12800
T=70000
((Ntot=Nstar+Nnomass))
#echo $Ntot
NAME=Nt${Ntot}_Ns${Nstar}_B${B}_Nn${Nnomass}_NBt_${T}
echo $NAME


rm *.info *.input *.10


./mcluster -N $Nstar -B $B -m 0.08 -m 10 -s 12345 -u $unit -C 5 -T $T -o $NAME > log.gen
var=`awk '$1=="scalingInfo" {print $4, $7}' log.gen`
echo $var
#python aei.py $var $NAME

echo 'first step finished!'

./mcluster -N $Nnomass -m 0.08 -m 10  -s 54321 -u $unit -C 5 -T $T  > test.log.gen
#./mcluster -N 128000 -n 128000 -B 12800 -m 0.08 -m 10 -s 54321 -t 3 -u 0 -C 5 -T 100 > test.log.gen
#exit
awk 'BEGIN {k = sqrt(0.5/0.5)} { print $1 * 0.000000001, $2, $3, $4, $5 * k, $6 * k, $7 * k}' test.dat.10 > dat.10
#awk 'BEGIN {k = sqrt(0.5/0.5)} { printf("%20.18f %20.18f %20.18f %20.18f %20.18f %20.18f %20.18f\n" $1 * 0.000000001, $2, $3, $4, $5 * k, $6 * k, $7 * k)}' test.dat.10 > dat.10
#rm ${NAME}COMs.*
rm test.*


echo 'second step finished!'

cat ${NAME}.dat.10 >> dat.10


rm -rf ../nbody && mkdir ../nbody
#mv dat.10 *.input ../nbody/
var=`awk 'NR==2,NR==2{print "'$Ntot'", $2,$3,$4,$5,$6,"'$Nnomass'"}' ${NAME}.input`
sed '2c '"$var"'' *.input > ../nbody/${NAME}.input
mv dat.10 ../nbody/
#cp dat.10 *.input start.lsf start.sh ../nbody/

