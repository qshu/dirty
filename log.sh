shopt -s extglob
workPath="dirty"
runPath="dirty"
logPath="LOG"
myHome="/home/qi/work"

cd ${myHome}/$workPath/ && echo 'enter home directory succesfully'
mkdir -p $logPath
#./configure --enable-mcmodel=large --with-par=b1m --enable-hdf5 > $logPath/log.configure
#./configure  --with-par=b32k --enable-hdf5 > $logPath/log.configure 
#./configure --enable-mcmodel=large  --with-par=b1m --enable-hdf5 > $logPath/log.configure 
#./configure --enable-mcmodel=large  --with-par=b1m > $logPath/log.configure 

#./configure  --with-par=b64k > $logPath/log.configure 
#NAME=nbody6++.avx.gpu.mpi

./configure  --with-par=b64k --disable-gpu > $logPath/log.configure 
NAME=nbody6++.avx.mpi

#make clean
rm -f ./build/*.o ./build/*.so ./build/nbody6++.avx* ./build/nb6++dumpb2a ./build/nb6++snapshot 
make -j > $logPath/log.make

str="rm !(*.input)";
cd $runPath/ && eval $str 
echo "********************** log.run begin ****************" > log.run
ls *.input >> log.run
echo "******* input file begin *******" >> log.run
cat *.input >> log.run
echo "******* input file end *********" >> log.run
echo "" >> log.run

echo "executable program: " $NAME >> log.run
date >> log.run
tstart=$(date +%s)
echo "******* simulation begin *******" >> log.run

# record the running log locally, keep logpath clean
#../build/nbody6++.avx.gpu.mpi < CometsProject.input >> log.run
../build/$NAME < CometsProject.input >> log.run

tend=$(date +%s)
echo "******* simulation end *********" >> log.run
date >> log.run
echo "" >> log.run

echo "Time cost in simulation: " $(($tend - $tstart)) "seconds" >> log.run
echo "********************** log.run end ******************" >> log.run

cp log.run ../$logPath/

