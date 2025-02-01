from generate_file import generate_solution, generate_file_name
from plot_file import load_solution
import os


solutions = []

d_ranges = {
    "cc": range(-95, 105, 10),
    "uc": range(-95, 105, 10),
    "ss": range(0, 90, 10),
    "ue": range(0, 60, 10),

}


def divide_by_100(n):
    return float(n)/100

def do_nothing(n):
    return n

d_process_ranges = {
    "cc": divide_by_100,
    "uc": divide_by_100,
    "ss": do_nothing,
    "rr": do_nothing,
    "ue": do_nothing,

}

for exp in d_ranges.keys():
    for gamma in [None, 1, 2, 3, 4, 5]:
        for lm in [None, 0.15]:
            for d in d_ranges[exp]:
                solutions.append([exp, gamma, lm, d])
            

for x in solutions:
    new_file = generate_file_name(x[0], x[1], x[2], d_process_ranges[x[0]](x[3]))
    print("\n\n" + new_file + "\n\n")

    if (os.path.isfile(new_file)):
        print(new_file + " already exists. Skipping.")

    else:
        print(new_file + " Doesn't exist yet.")
        solution = generate_solution(x[0], x[0], x[0], 
                                     Gamma=x[1], lambd=x[2], 
                                     timesteps=50,
                                     deformation=d_process_ranges[x[0]](x[3]))
        print("\n\n" + solution + " GENERATED\n\n")
        load_solution(solution)

print("All listed files exist!")


