# make dir
mkdir -p ./data
mkdir -p ./figure
mkdir -p ./bin

python gap_eval.py plot

# compile and run lp_eval
g++ -O2 ./cpp_eval/main.cpp -o ./bin/lp_eval
./bin/lp_eval minimum_weight 6 4 5 100000
./bin/lp_eval minimum_weight 6 4 7 100000
./bin/lp_eval minimum_weight 6 8 5 100000
./bin/lp_eval minimum_weight 6 8 7 100000
# ./bin/lp_eval minimum_weight 10 4 5 10000
# ./bin/lp_eval minimum_weight 10 4 7 10000
# ./bin/lp_eval minimum_weight 10 8 5 10000
# ./bin/lp_eval minimum_weight 10 8 7 10000
