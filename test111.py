import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
#%%
img=Image.open('SimWakeInDas.png')
print(img.size)
#img=np.array(img)

plt.imshow(img)
plt.show()
#%%