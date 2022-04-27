setenv JULIA_DEPOT_PATH /usr/local/usrapps/gluonsaturation/julia
path=`pwd`

for L in 16 24 32 48 ; do

 id=$RANDOM # not important -- just to provide a random seed for the start

 # create a tmp file for submission 
 TMPFILE=`mktemp tmp.XXXXXXXXXXXX`
 
 # populate teh file with needed script 
 cp run_short.sh $TMPFILE
 echo "julia -t 16 modelA.jl  $id $L $RANDOM  >  /rsstu/users/v/vskokov/gluon/criticaldynamic/tmp/rew_${L}_${ser}_${id}.dat"  >> $TMPFILE
 echo "rm $path/$TMPFILE "  >> $TMPFILE
 
 # submit
 cp run_short.sh $TMPFILE
 chmod u+x $TMPFILE
 bsub < $TMPFILE

done
