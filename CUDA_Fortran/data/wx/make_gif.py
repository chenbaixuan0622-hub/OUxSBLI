from PIL import Image
import glob

files = sorted(glob.glob('./*.png'))
images = list(map(lambda file : Image.open(file), files))
images[0].save('Vorticity.gif', save_all = True, append_images = images[1:], duration = 10, loop = 0)
