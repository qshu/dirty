shopt -s extglob
workPath="dirty"
logPath="LOG"


mkdir -p $logPath
#./configure --enable-mcmodel=large --with-par=b1m --enable-hdf5 > $logPath/log.configure
./configure  --with-par=b32k --enable-hdf5 > $logPath/log.configure 
#make clean
rm -f ./build/*.o ./build/*.so ./build/nbody6++.avx* ./build/nb6++dumpb2a ./build/nb6++snapshot 
make -j > $logPath/log.make

str="rm !(*.input)";
cd $workPath/ && eval $str 
echo "********************** log.run begin ****************" > log.run
ls *.input >> log.run
echo "******* input file begin *******" >> log.run
cat *.input >> log.run
echo "******* input file end *********" >> log.run
echo "" >> log.run

tstart=$(date +%s)
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
echo $(($tend - $tstart)) "seconds" >> log.run
echo "********************** log.run end ******************" >> log.run

mv log.run ../$logPath/

