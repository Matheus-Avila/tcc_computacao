import os
import moviepy.video.io.ImageSequenceClip
image_folder='figs'
fps=0.5

#image_files = [image_folder+'/'+img for img in os.listdir(image_folder) if img.endswith(".png")]
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(['0.0000.png','1.0000.png','2.0000.png','3.0000.png','4.0000.png','5.0000.png','6.0000.png','6.9990.png','0.0000.png','0.0000.png',], fps=fps)
clip.write_videofile('my_video.mp4')