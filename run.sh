# make dir
mkdir -p ./data
mkdir -p ./figure
mkdir -p ./bin

# run gap_eval.py
# python gap_eval.py run 4 4 3 10000
python gap_eval.py run 4 4 5 100000
python gap_eval.py run 4 4 7 100000
# python gap_eval.py run 4 4 9 10000
python gap_eval.py plot

# python gap_eval.py run 1 10 3 10000
python gap_eval.py run 1 10 5 100000
python gap_eval.py run 1 10 7 100000
# python gap_eval.py run 1 10 9 10000
python gap_eval.py plot

# compile and run lp_eval
g++ -O2 ./cpp_eval/main.cpp -o ./bin/lp_eval
# ./bin/lp_eval minimum_weight 6 4 3 1000
./bin/lp_eval minimum_weight 6 4 5 1000
./bin/lp_eval minimum_weight 6 4 7 1000
# ./bin/lp_eval minimum_weight 6 4 9 1000
