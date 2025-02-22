from lagrange_gen import Experiment, TIMESTEPS
import threading


IMAGE_COUNT = 5
IMAGES = range(0, TIMESTEPS, round(TIMESTEPS/IMAGE_COUNT))

TEMPS = []
for x in range(-30, 2, 2):
    TEMPS.append(x)

EXP_TYPES = [ "cc", "ue", "uc", "ss"]
LAMS = [False, True]



def run(action, exptypes = EXP_TYPES, temps = TEMPS):
    for exp in exptypes:
        for tem in temps:
            print(exp, tem)
            
            if exp != "ss" :
                e1 = Experiment(exp, "zz", temp = tem) 
            else:
                e1 = Experiment(exp, "xz", temp = tem) 

            print(f"GAMMA: {e1.Gamma}, LAMD: 0.15, TEMP: {e1.temp}, EXP: {e1.exptype}")
            action(e1)


def file_insides(experiment):

    print(experiment.get_dataframe())


def main():
    run(file_insides, ['ss'], [-30])

main()