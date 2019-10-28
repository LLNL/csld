#!/bin/bash
set -e
if [ $# -lt 11 ]; then
    echo "Usage: $0 prim_poscar potcar prim_kpoints prim_incar sc.txt sc_kpoints sc_incar \"scale1(relative to prim_poscar) [scale2...]\" disp1 ndisp1 spec_disp1 [disp2 ndisp2 spec_disp2 ...]"
    echo "suggetion: mkkp.sh prim_poscar|sc_poscar 3000 > prim_kpoints|sc_kpoints"
    echo "Example:"
    echo '../prepare-volume.sh POSCAR  POTCAR  KPOINTS_prim INCAR_prim  ../sc.txt  KPOINTS INCAR "4.01 4.08 4.13"   0.01 2 "4 5.5" 0.05 1 "5 5.3" 0.12 3 "6 5.2"  0.3 2 "6 5.1"'
    exit 0
fi
prim_poscar=$1
shift
prim_kpoints=$1
shift
potcar=$1
shift
prim_incar=$1
shift
sc=$1
shift
sc_kpoints=$1
shift
sc_incar=$1
shift
scale_list=$1
shift

#echo "debug prim_poscar prim_kpoints prim_incar sc.txt sc_kpoints sc_incar scales $prim_poscar $prim_kpoints $prim_incar $sc $sc_kpoints $sc_incar $scale_list"
#exit 0

a0=`echo $scale_list |awk '{print $1}'`
run0=a$a0
for a in $scale_list; do
    run=a$a
    echo "Generating displacements for a =  $a"
    prim_dir=$run/prim
    mkdir -p $run
    mkdir -p $prim_dir
    awk -v a=$a '{if (NR==2) $1=a; print}' $prim_poscar > $prim_dir/POSCAR
    ln -sr $prim_kpoints $prim_dir/KPOINTS
    ln -sr $prim_incar $prim_dir/INCAR
    ln -sr $potcar  $prim_dir/POTCAR
    polaron_main --task supercell --p1 $sc --prim $prim_dir/POSCAR > $run/SPOSCAR
    if [ $a == $a0 ]; then
	# create displacements
	while (( "$#" )); do
	    disp=$1; shift; ndisp=$1; shift; spec=$1; shift
	    polaron_main --task rand_disp_dir -R $disp --npolaron $ndisp --p1 $run/SPOSCAR --p2 "$spec"
	    mv dir_*/ $run
	done
	for i in $run/dir_*; do
	    polaron_main --task cmp_pos --p1 $run/SPOSCAR  --p2 $i/POSCAR  --summary > $i/disp.txt
	done
    else
	for i in $run0/dir_*; do 
	    #	    polaron_main --task vol_disp --p1 POSCAR  --p2 $i/POSCAR  -R `python -c "print($scale/$a0)"`  > $dir/$i/POSCAR
	    dir=$run/${i#*/}
	    mkdir -p $dir
	    polaron_main --task add_disp --p1 $run/SPOSCAR --p2 $i/disp.txt > $dir/POSCAR
	done
    fi
    ln -sr $sc_kpoints $dir/KPOINTS
    ln -sr $sc_incar $dir/INCAR
    ln -sr $potcar $dir/POTCAR
done
