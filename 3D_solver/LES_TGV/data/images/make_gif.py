from PIL import Image
import glob

files = sorted(glob.glob('./*.png'))  
images = list(map(lambda file : Image.open(file) , files))
images[0].save('Q0.1.gif' , save_all = True , append_images = images[1:] , duration = 50 , loop = 0)

