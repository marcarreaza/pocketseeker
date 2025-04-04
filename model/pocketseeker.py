import sys
sys.path.append('../extract_features')

from concativity.concavity_feature import concavity_feature
from core_distance.core_distance_feature import core_distance
from physicochemical.physicochemical_features import physicochemical_feature
from PSSM.PSSM_feature import PSSM_feature

a = "protein.pdb"
concavity = concavity_feature(a)
distance_to_core = core_distance(a)
physicochemical = physicochemical_feature (a)
PSSM = PSSM_feature (a)
