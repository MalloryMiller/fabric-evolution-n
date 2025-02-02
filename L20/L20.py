from experiment import Experiment 

for x in [14, 16, 18, 20]:

    e1 = Experiment("cc", "zz", temp = -30, truncation=x) 

    e1.generate_file()
    e1.load_solution([49])