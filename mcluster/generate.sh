make clean 
make mcluster
#./mcluster -N 10000 -B 1000 -s 12345 -u 1 -C 5 -o CometsProject



unit=0 # 0 --> Nbody unit;  1 --> astrophysics unit
Nstar=500000 # Ns: total number of particles or stars(for ffc)
B=0.1
Nnomass=50000
T=10
t=1 #standard solar neighbour tidal field
EPOCH=0  # Evolutionary epoch, in unit of Myr
((Ntot=Nstar+Nnomass))
#echo $Ntot
m=100

#orbital parameter:
amin=10
amax=1000

NAME=Nt${Ntot}_Ns${Nstar}_B${B}_Nn${Nnomass}_NBt_${T}
echo $NAME
echo $Nnomass


rm *.info *.input *.10


./mcluster -N $Nstar -m 0.08 -m $m -s 1000 -u $unit -C 5 -T $T -t $t -o $NAME > log.gen

var=`awk '$1=="scalingInfo" {print $4, $7}' log.gen`
echo $var
#python aei.py $var $NAME

echo 'first step finished!'

echo $var $NAME $Nstar $Nnomass
python aei.py $var $NAME $Nstar $Nnomass $amin $amax

#./mcluster -N $Nnomass -B 0 -m 0.08 -m $m  -s 2001 -u $unit -C 5 -nbin 2 -T $T -t $t > test.log.gen
#./mcluster -N 128000 -n 128000 -B 12800 -m 0.08 -m 10 -s 54321 -t 3 -u 0 -C 5 -T 100 > test.log.gen
#exit
#awk 'BEGIN {k = sqrt(0.5/0.5)} { print $1 * 0.000000001, $2, $3, $4, $5 * k, $6 * k, $7 * k}' test.dat.10 >> ${NAME}.dat.10
#awk 'BEGIN {k = sqrt(0.5/0.5)} { printf("%20.18f %20.18f %20.18f %20.18f %20.18f %20.18f %20.18f\n" $1 * 0.000000001, $2, $3, $4, $5 * k, $6 * k, $7 * k)}' test.dat.10 > dat.10
#rm ${NAME}COMs.*
#rm test.*

#mv test.dat.10 ${NAME}.dat.10

echo 'second step finished!'

#mv ${NAME}.dat.10  dat.10

rm -rf ../nbody && mkdir ../nbody
#cp dat.10 *.input ../nbody/
var1=`awk 'NR==2,NR==2{print "'$Ntot'", $2,$3,$4,$5,$6,"'$Nnomass'"}' ${NAME}.input`
sed '2c '"$var1"'' ${NAME}.input > ../nbody/${NAME}.input
mv ../nbody/${NAME}.input ./
var2=`awk 'NR==10,NR==10{print $1,$2,$3,"'$Nnomass'", $5,$6,$7,$8}' ${NAME}.input`
sed '10c '"$var2"'' ${NAME}.input > ../nbody/${NAME}.input
cp dat.10 ../nbody/
#cp dat.10 *.input start.lsf start.sh ../nbody/

