for fp_len in 64 128 256 512 1024 ; do
    echo ${fp_len}
    sed 's/@fp_len@/'"${fp_len}"'/g'   main_1_model_train_model.py  >  main_1_model_train.py
    python main_1_model_train.py  &>  OUTPUT_${fp_len}.txt
done
