from experiment import Experiment 


TIMESTEPS = 1501

IMAGE_COUNT = 10
IMAGES = range(0, TIMESTEPS, round(TIMESTEPS/IMAGE_COUNT))

TEMPS = [False]
for x in range(-30, 2, 2):
    TEMPS.append(x)

EXP_TYPES = [ "cc", "ue", "uc", "ss"]
LAMS = [False, True]
L = 20



def main(action):
    for exp in EXP_TYPES:
        for tem in TEMPS:
            for la in LAMS:
                
                if exp != "ss" :
                    e1 = Experiment(exp, "zz", temp = tem, lambd=la, truncation=L, timesteps=TIMESTEPS) 
                else:
                    e1 = Experiment(exp, "xz", temp = tem, lambd=la, truncation=L, timesteps=TIMESTEPS) 

                print(f"GAMMA: {e1.Gamma}, LAMD: {e1.lamd}, TEMP: {e1.temp}, EXP: {e1.exptype}")
                action(e1)


def create_ncs(experiment):
    experiment.generate_file(remake=False)

def create_images(experiment):
    experiment.load_solution(IMAGES)


#main(create_ncs)
main(create_images)