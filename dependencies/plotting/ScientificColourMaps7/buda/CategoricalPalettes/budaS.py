# 
#         budaS
#                   www.fabiocrameri.ch/colourmaps
from matplotlib.colors import LinearSegmentedColormap      
      
cm_data = [[0.70015, 0.0027445, 0.70061],      
           [1, 1, 0.4002],      
           [0.80548, 0.52274, 0.495],      
           [0.73557, 0.30688, 0.56485],      
           [0.86025, 0.74, 0.44036],      
           [0.7723, 0.41654, 0.52814],      
           [0.89165, 0.85727, 0.41182],      
           [0.70312, 0.18374, 0.61179],      
           [0.83232, 0.62832, 0.46822],      
           [0.87466, 0.79754, 0.42602],      
           [0.93038, 0.92449, 0.40233],      
           [0.78971, 0.46941, 0.51074],      
           [0.71634, 0.24759, 0.58531],      
           [0.70042, 0.11127, 0.65003],      
           [0.84616, 0.68368, 0.45442],      
           [0.75431, 0.36274, 0.54612],      
           [0.81861, 0.57361, 0.4819],      
           [0.70095, 0.14919, 0.62928],      
           [0.76337, 0.38982, 0.53706],      
           [0.83921, 0.6559, 0.46134],      
           [0.72586, 0.27775, 0.57472],      
           [0.85317, 0.7117, 0.44741],      
           [0.88228, 0.82696, 0.41878],      
           [0.79789, 0.49592, 0.50257],      
           [0.8674, 0.76859, 0.43323],      
           [0.9064, 0.88953, 0.40601],      
           [0.70826, 0.21633, 0.59734],      
           [0.82545, 0.6009, 0.47508],      
           [0.70034, 0.066138, 0.67402],      
           [0.96301, 0.96169, 0.40079],      
           [0.74506, 0.33517, 0.55536],      
           [0.78111, 0.44302, 0.51934],      
           [0.81258, 0.54981, 0.48792],      
           [0.80176, 0.5093, 0.49871],      
           [0.98129, 0.98078, 0.40045],      
           [0.71201, 0.2321, 0.59107],      
           [0.70173, 0.1668, 0.62015],      
           [0.75886, 0.37635, 0.54156],      
           [0.78545, 0.4562, 0.51501],      
           [0.79388, 0.48264, 0.50659],      
           [0.7497, 0.34902, 0.55071],      
           [0.8567, 0.72581, 0.44391],      
           [0.70057, 0.13079, 0.63925],      
           [0.80908, 0.53625, 0.49142],      
           [0.74035, 0.32111, 0.56008],      
           [0.76785, 0.40323, 0.53259],      
           [0.70037, 0.090118, 0.66165],      
           [0.86382, 0.75426, 0.43681],      
           [0.77673, 0.4298, 0.52372],      
           [0.72101, 0.2628, 0.5799],      
           [0.70028, 0.036129, 0.68707],      
           [0.82888, 0.61459, 0.47165],      
           [0.87101, 0.78302, 0.42964],      
           [0.91714, 0.90667, 0.40388],      
           [0.70527, 0.2002, 0.60421],      
           [0.94585, 0.94289, 0.40136],      
           [0.82204, 0.58724, 0.47848],      
           [0.88657, 0.84195, 0.41522],      
           [0.87838, 0.81218, 0.42239],      
           [0.89805, 0.87309, 0.40869],      
           [0.73073, 0.29243, 0.56973],      
           [0.84268, 0.66976, 0.45788],      
           [0.84966, 0.69765, 0.45092],      
           [0.83576, 0.64209, 0.46479],      
           [0.81518, 0.56, 0.48532],      
           [0.99068, 0.9904, 0.40032],      
           [0.79985, 0.5026, 0.50062],      
           [0.84791, 0.69065, 0.45268],      
           [0.77451, 0.42317, 0.52593],      
           [0.70035, 0.078557, 0.66775],      
           [0.74271, 0.32815, 0.55771],      
           [0.87283, 0.79027, 0.42783],      
           [0.76112, 0.3831, 0.53931],      
           [0.70073, 0.14008, 0.63417],      
           [0.81084, 0.54303, 0.48964],      
           [0.77008, 0.4099, 0.53036],      
           [0.86203, 0.74711, 0.43858],      
           [0.8169, 0.56681, 0.4836],      
           [0.8803, 0.81954, 0.42058],      
           [0.88436, 0.83442, 0.41699],      
           [0.71411, 0.23986, 0.58813],      
           [0.70234, 0.17534, 0.61587],      
           [0.8765, 0.80485, 0.4242],      
           [0.83404, 0.6352, 0.46651],      
           [0.70023, 0.019196, 0.69378],      
           [0.91146, 0.89801, 0.40486],      
           [0.80364, 0.516, 0.49683],      
           [0.82032, 0.58043, 0.48019],      
           [0.8306, 0.62145, 0.46994],      
           [0.85847, 0.73289, 0.44214],      
           [0.89465, 0.86511, 0.41022],      
           [0.75659, 0.36956, 0.54384],      
           [0.70128, 0.15808, 0.62461],      
           [0.85141, 0.70467, 0.44918],      
           [0.74739, 0.34211, 0.55303],      
           [0.88898, 0.84956, 0.4135],      
           [0.75202, 0.35591, 0.54842],      
           [0.86921, 0.7758, 0.43144],      
           [0.71005, 0.22422, 0.59413],      
           [0.70409, 0.19203, 0.60791]]      
      
budaS_map = LinearSegmentedColormap.from_list('budaS', cm_data)      
# For use of "viscm view"      
test_cm = budaS_map      
      
if __name__ == "__main__":      
    import matplotlib.pyplot as plt      
    import numpy as np      
      
    try:      
        from viscm import viscm      
        viscm(budaS_map)      
    except ImportError:      
        print("viscm not found, falling back on simple display")      
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',      
                   cmap=budaS_map)      
    plt.show()      
