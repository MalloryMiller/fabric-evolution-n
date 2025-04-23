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
    '''
    creates a GeneratedExperiment for every experiment combination of the given
    exptypes and temps. Then runs the given action() function with the generated exeriment as
    the only input to the function.
    '''
    for exp in exptypes:
        for tem in temps:
            print(exp, tem)
            
            if exp != "ss" :
                e1 = GeneratedExperiment(exp, "zz", temp = tem) 
            else:
                e1 = GeneratedExperiment(exp, "xz", temp = tem) 

            print(f"GAMMA: {e1.Gamma}, LAMD: 0.15, TEMP: {e1.temp}, EXP: {e1.exptype}")
            action(e1)


def create_ncs(experiment, remake=False):
    '''
    Creates nc files for the experiment input. 
    By default it does not overwrite preexisting nc files,
    and if one already exists the function will return after printing a warning.
    '''
    experiment.generate_file(remake)


def create_images(experiment):
    '''
    Creates an image files for the experiment input. 
    An nc file must be generated for this experiment already for this function to work.
    By default it does not overwrite a preexisting image folder,
    and if one already exists the function will return after printing a warning.
    '''
    experiment.load_solution([TIMESTEPS], isolate=True, remake = True)


def main():
    '''
    Creates nc files for all data in the experimental range
    and saves images corresponding to those experiments.
    '''
    run(create_ncs, EXP_TYPES, TEMPS)
    run(create_images, EXP_TYPES, TEMPS)

main()