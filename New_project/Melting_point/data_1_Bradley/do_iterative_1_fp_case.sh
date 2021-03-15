for fp_case in $(seq 0 7) ; do
    echo ${fp_case}
    sed 's/@fp_case@/'"${fp_case}"'/g'   main_1_model_train_model.py  >  main_1_model_train.py
    python main_1_model_train.py  &>  OUTPUT_${fp_case}.txt
done
