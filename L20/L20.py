from experiment_generation import GeneratedExperiment, TIMESTEPS
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
                e1 = GeneratedExperiment(exp, "zz", temp = tem) 
            else:
                e1 = GeneratedExperiment(exp, "xz", temp = tem) 

            print(f"GAMMA: {e1.Gamma}, LAMD: 0.15, TEMP: {e1.temp}, EXP: {e1.exptype}")
            action(e1)


def create_ncs(experiment):
    experiment.generate_file(remake=False)


def create_images(experiment):
    experiment.load_solution([TIMESTEPS], isolate=True, remake = True)


def main():
    run(create_ncs, EXP_TYPES, TEMPS)
    #run(create_images, EXP_TYPES, TEMPS)

main()