TASKNUM=4
x=1
while [ $x -le $TASKNUM ]
do
  ./isosurface gs.bp iso.bp 0.5
  x=$(( $x + 1 ))
done