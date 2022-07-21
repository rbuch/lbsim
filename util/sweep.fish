#! /usr/bin/env fish

set pmax $argv[1]
set omax $argv[2]
set smax $argv[3]

set dir "run-"$pmax"-"$omax"-"$smax"_"(date +%s)
mkdir $dir

for i in (seq $argv[1])
    for j in (seq $argv[2])
	set procs (math 2^$i)
	set objs (math 10^$j)

	if $objs <= $procs
	    continue
	end
	
	set conf $procs"_"$objs
	echo $conf
	mkdir $dir/$conf
	for s in (seq 0 4 $argv[3])
	    for k in (seq 0 3)
		set seed (math $s + $k)
		set logfile $procs"_"$objs"_"$seed".log"
		#echo "./lbsim (math 2^$i) (math 10^$j) $seed & > $dir/$conf/$logfile"
		./lbsim $procs $objs $seed > $dir/$conf/$logfile &
	    end
	    wait
	end
    end
end
