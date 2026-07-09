# make dir
mkdir -p ./data
mkdir -p ./figure
mkdir -p ./bin

# run gap_eval.py
python gap_eval.py run 4 4 5 100000
python gap_eval.py run 4 4 7 100000

python gap_eval.py run 4 8 5 100000
python gap_eval.py run 4 8 7 100000

python gap_eval.py run 8 4 5 100000
python gap_eval.py run 8 4 7 100000

python gap_eval.py run 8 8 5 100000
python gap_eval.py run 8 8 7 100000
