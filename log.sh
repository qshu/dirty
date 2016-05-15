shopt -s extglob
workPath="dirty"
logPath="LOG"

mkdir -p $logPath
./configure --enable-hdf5 > $logPath/log.configure
make clean
make -j 16 > $logPath/log.make

str="rm !(*.input)";
cd $workPath/ && eval $str 
echo "********************** log.run begin ****************" > log.run
ls *.input >> log.run
echo "******* input file begin *******" >> log.run
cat *.input >> log.run
echo "******* input file end *********" >> log.run
echo "" >> log.run

tstsrt=$(date +%s)
echo "******* simulation begin *******" >> log.run
date >> log.run
echo "********************************" >> log.run

../build/nbody6++.avx* < *.input >> log.run

tend=$(date +%s)
echo "******* simulation end *********" >> log.run
date >> log.run
echo "********************************" >> log.run

echo "" >> log.run
echo "Time cost in simulation:" >> log.run
echo $(($tend - $tstsrt)) "seconds" >> log.run
echo "********************** log.run end ******************" >> log.run

mv log.run ../$logPath/

