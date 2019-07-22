# Edit the number of parallel executions. 
N=6

g++  -std=c++11 src/*.cpp src/*.h -o binary -O3
echo "compiled"


open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    printf '%.3d' $? >&3
    )&
}


run_one(){
    declare -i TIME=$SECONDS
    INSTANCE_NAME=$(basename "$f" .dat)
    SEED=$1    
    POPSIZE=972
    EXPONENT=5.14
    f=$2
    sep="_"
    #FILE="out/$INSTANCE_NAME$sep$SEED.txt"
    FILE="results.txt"    
    echo "Working on $INSTANCE_NAME"
    ./binary $INSTANCE_NAME $SEED $POPSIZE $EXPONENT < $f >> $FILE
    echo "" >> $FILE
    # expr $SECONDS - $TIME >> $FILE ;
    # sleep $SEED
}


open_sem $N
for f in QAPinstances/*.dat; do
for SEED in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21; do
    run_with_lock run_one $SEED $f;
done 
done

