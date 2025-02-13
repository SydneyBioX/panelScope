import copy, os
import random
from collections import deque
import numpy as np
import pandas as pd

class Gene_panel:
    """
    Parameters:
    config : Gene set, str  
    info : Save information of Gene_panel, str
    result : The metric of a Gene_panel, float
    """
    def __init__(self, config=None, info=None, result=None):
        self.config = config
        self.result = result
    def __str__(self):
        return "config :" + str(self.config) + ", result: " + str(self.result)


def random_config_generator(search_space, is_rand, random_state, old_config=None):
    """
    Parameters:
    search_space: dict like {gene_position_idx:[all potential gene list]
    is_rand: dict of which position should be randomly pick new gene
    random_state: a random seed
    old_config: a old config
    """
    if old_config:
        new_config=copy.deepcopy(old_config)
    else:
        new_config={item:None for item in search_space}  
    for item in is_rand:
        if is_rand[item]:
            _index = random_state.randint(len(search_space[item]))
            new_config[item] = search_space[item][_index]
    return new_config

class Evolution:
    """
    Parameters:
    optimize_mode: the direction of optimization 'maximize' or 'minimize', str
    population_size: The initial size of the evolution population, int
    """

    def __init__(self, panel_score_path, optimize_mode='maximize', population_size=32, objmode='overall'):
        self.panel_score_path = panel_score_path
        self.optimize_mode = optimize_mode
        self.population_size = population_size
        self.objmode = objmode

        self.running_trials = {}
        self.num_running_trials = 0
        self.random_state = None
        self.population = None
        self.space = None
        self.unsatisfied = 0
        self.param_ids = deque()

    def update_search_space(self, search_space):
        """
        Parameters:
        search_space : dict like {gene_position_idx:[all potential gene list]}
        """
        self.space = search_space
        self.random_state = np.random.RandomState()
        self.population = []

        for _ in range(self.population_size):
            self._random_generate_Gene_panel()   
    def _random_generate_Gene_panel(self):
        is_rand = dict()
        for item in self.space:
            is_rand[item] = True
        config = random_config_generator(self.space, is_rand, self.random_state)
        self.population.append(Gene_panel(config=config))

    def generate_multiple_panels(self, parameter_id_list, **kwargs):
        """
        Parameters:
        parameter_id_list : Unique index for each panel, list of int
        
        Returns:
        A list of newly generated configurations
        """
        result = []
        for parameter_id in parameter_id_list:
            res = self.generate_panel(parameter_id, **kwargs)
            self.num_running_trials += 1
            result.append(res)
        return result
    
    
    def generate_panel(self, parameter_id, **kwargs):
        """
        function for generate_multiple_panels
        Parameters:
        parameter_id : int
        """
        if not self.population:
            raise RuntimeError('The population is empty')

        if self.num_running_trials >= self.population_size:
            print("Population_size is suggested to be larger than Concurrency")
            self.unsatisfied += 1
            self.param_ids.append(parameter_id)
            raise RuntimeError('Population_size is suggested to be larger than Concurrency')

        return self._generate_Gene_panel(parameter_id)  
    def _generate_Gene_panel(self, parameter_id):
        pos = -1
        for i in range(len(self.population)):
            if self.population[i].result is None:
                pos = i
                break
        if pos != -1:
            indiv = copy.deepcopy(self.population[pos])
            self.population.pop(pos)
        else:
            random.shuffle(self.population)
            if len(self.population) > 1 and self.population[0].result < self.population[1].result:
                self.population[0] = self.population[1]
            space = list(self.population[0].config.keys())
            is_rand = dict()
            mutation_pos = space[random.randint(0, len(space)-1)]
            for i in self.space.keys():
                is_rand[i] = (i == mutation_pos)
            config = random_config_generator(
                self.space, is_rand, self.random_state, self.population[0].config)
            if len(self.population) > 1:
                self.population.pop(1)
            indiv = Gene_panel(config=config)
        self.running_trials[parameter_id] = indiv
        return indiv.config
    
    def receive_trial_result(self, parameter_id, **kwargs):
        # try:
        scores = self._extract_scalar_reward(parameter_id)      
        if self.objmode=="overall":
            reward=scores['overall']
        elif self.objmode=='cts':
            reward=scores['cts']
        elif self.objmode=='pathway':
            reward=scores['pathway']
        elif self.objmode=='spatial':
            reward=scores['spatial']
        elif self.objmode=='tv':
            reward=scores['tv']
        else:
            reward=scores['corr']
                    
        if parameter_id not in self.running_trials:
            raise RuntimeError('Received parameter_id %s not in running_trials.', parameter_id)
        config = self.running_trials[parameter_id].config
        self.running_trials.pop(parameter_id)
        if self.optimize_mode != 'maximize':
            reward = -reward
        indiv = Gene_panel(config=config, result=reward)
        self.population.append(indiv)
        with open(self.objmode+'_reward.txt','a') as ffff:
            ffff.write(str(reward)+'\n')
        return [True, reward]
        # except:
        #     return [False, None]
            
    def trial_end(self, parameter_id, success, **kwargs):
        self.num_running_trials -= 1
        if not success:
            self.running_trials.pop(parameter_id)
            self._random_generate_Gene_panel()
        if self.unsatisfied > 1:
            param_id = self.param_ids.popleft()
            config = self._generate_Gene_panel(param_id)
            self.unsatisfied -= 1
            self.num_running_trials += 1

    def _extract_scalar_reward(self, param_id):
        scores_df = pd.read_csv(os.path.join(self.panel_score_path, str(param_id)+'_scores.csv'))
        scores={}
        scores['cts'] = scores_df["Celltype_specificity_score"].iloc[0]
        scores['corr'] = scores_df["Correlation_ratio"].iloc[0]
        scores['pathway'] = scores_df["Pathway_score"].iloc[0]
        scores['spatial'] = scores_df["Spatial_score"].iloc[0]
        scores['tv'] = scores_df["Transcriptional_variability"].iloc[0]
        scores['overall'] = scores['cts'] + scores['corr'] + scores['pathway'] + scores['spatial'] + scores['tv']
        return scores