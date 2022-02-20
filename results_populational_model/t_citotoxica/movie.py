import os
import moviepy.video.io.ImageSequenceClip
image_folder='figs'
fps=0.5

#image_files = [image_folder+'/'+img for img in os.listdir(image_folder) if img.endswith(".png")]
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(['fig0.png','fig0.1000.png','fig0.1000.png','fig0.2000.png','fig0.3000.png','fig0.4000.png','fig0.4990.png','fig0.4990.png','fig0.4990.png',], fps=fps)
clip.write_videofile('my_video.mp4')