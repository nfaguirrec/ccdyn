#!/bin/bash 

IFILE=$1
NTRAJ=$2

IFILE_EFF="$IFILE-eff$$"

runDroplet(){
	local idLabel=`echo $1 | awk '{printf "%04d\n",$1}'`
	
	sed 's/\&data/\&data\nidLabel='$idLabel'/g' $IFILE_EFF > $idLabel.in
	HeDropletDynamics $idLabel.in > $idLabel.out
}

main(){
	local i=0
	local j=0
	local ij=0
	local nProcs=1
	local scratch=""
	local work=""
	
	if [ -n "$PBS_NODEFILE" ]
	then
		nProcs=`cat "$PBS_NODEFILE" 2> /dev/null | wc -l`
		scratch=$pbstmp
		work=$PBS_O_WORKDIR
	else
#		nProcs=`cat /proc/cpuinfo | grep processor | wc -l`
		nProcs=2
		scratch=/scratch/$USER/HDD/$IFILE_EFF
		work=$PWD
		mkdir -p $scratch
	fi
	
	cd $work
	
	dep=`awk 'BEGIN{FS="="}($0~"initialDensity[[:blank:]]*="){print $2}' $IFILE | sed -r 's/[[:blank:]\"]+//g'`
	dep="$PWD/$dep"
	awk '{ if($0~"initialDensity[[:blank:]]*=") print "initialDensity=\"'$dep'\""; else print $0 }' $IFILE > $IFILE_EFF
	
	mv $IFILE_EFF $scratch
	
	cd $scratch
	
	for i in `seq 0 $(( $NTRAJ/$nProcs-1 ))`
	do
		echo -n "Running trajectories ` seq -s " " $(( 1+$nProcs*$i )) $(( $nProcs+$nProcs*$i )) ` ... "
		
		for j in `seq 1 $nProcs`
		do
			ij=$(( $j+$nProcs*$i ))
			runDroplet $ij &
		done
		
		wait
		echo "OK"
	done
	
	rm $IFILE_EFF
	cp -r * $work/
	rm -rf $scratch
}

main $*
