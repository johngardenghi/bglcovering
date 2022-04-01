#!/bin/bash

problem=$1
nvert=$2

if [ -z "$problem" ]; then
    echo "Usage: $0 america|cesaro|star|stoyan|regpol [nvert]"
    exit 1
fi

make

backupdir=~/documents/covering/tests/
balls=$(seq 10 10 100)
budget=(0 100 200 500 1000)

if [ ! -d "$backupdir" ]; then
    mkdir -p $backupdir
fi

for ((i=1; i<${#budget[@]}; i++)); do
    for nballs in ${balls[@]}; do
	cp ${backupdir}/bestrad-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$((i-1))]}).txt bestrad.txt

	echo "Running nballs = $nballs, itrial from $((${budget[$((i-1))]}+1)) to ${budget[$i]}, saving to ${backupdir}"

	for itrial in $(seq $((${budget[$((i-1))]}+1)) 1 ${budget[$i]}); do
	    if [ "$problem" = "regpol" ]; then
    		printf "regpol\n${nballs}\n${itrial}\nT\nF\n0.03\n0.15\n${nvert}\n" | timeout --preserve-status 40 ${PWD}/covering > /dev/null
		# printf "regpol\n${nballs}\n${itrial}\nT\nT\n${nvert}\n" | timeout --preserve-status 40 ${PWD}/covering > /dev/null
	    elif [ "$problem" = "america" ] || [ "$problem" = "cesaro" ] || [ "$problem" = "star" ] || [ "$problem" = "stoyan" ]; then
		printf "${problem}\n${nballs}\n${itrial}\nT\nF\n0.03\n0.15\n" | timeout --preserve-status 160 ${PWD}/covering > /dev/null
	    else
		echo "Usage: $0 america|cesaro|star|stoyan|regpol [nvert]"
		exit 1
	    fi

	    # Save iteration information
	    export statdtris=$?
	    if [ -e tabline.txt ] ; then
		sed -e 's/$/  BR/' -i tabline.txt
		cat tabline.txt >> table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt
		bestrad=$(cat bestrad.txt)
		bestitrial=$itrial
	    else
		if [ $statdtris -eq 0 ] && [ -e alsolver-interrupted-tabline.txt ]; then
		    sed -e 's/$/  S/' -i alsolver-interrupted-tabline.txt
		    cat alsolver-interrupted-tabline.txt >> table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt ;
		elif [ $statdtris != 0 ] && [ -e alsolver-interrupted-tabline.txt ]; then
		    sed -e 's/$/  F/' -i alsolver-interrupted-tabline.txt
		    cat alsolver-interrupted-tabline.txt >> table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt ;
		else
		    printf " no information\n" >> table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt ;
		fi
	    fi
	    rm -f tabline.txt alsolver-interrupted-tabline.txt \
	       solver-interrupted-tabline.txt newtkkt-interrupted-tabline.txt \
	       picture-solution-*
	done

	mv bestrad.txt ${backupdir}/bestrad-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt

	if [ "$backupdir" != "." ]; then
	    if [ "${budget[$i]}" = "0" ]; then
		mv table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt ${backupdir}
	    else
		cat table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt \
		    >> ${backupdir}/table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt

		rm -f table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt
	    fi
	fi

	# read -p "Press any key to resume ..."
    done
done
