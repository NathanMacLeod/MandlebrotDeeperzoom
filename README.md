More sophisticated mandlebrot program I wrote around May 2024. Compared to my previous mandlebrot program, this one features a software float implementation that allows for deeper zooms. 

Rather than using the software float for all pixels, I used a technique that allows calculating the orbit of a few points per frame, and then the orbits of the neighboring pixels is approximated. That method is described here:

https://dirkwhoffmann.github.io/DeepDrill/docs/Theory/Perturbation.html

I also use another technique for smoother coloring, described here:

https://dirkwhoffmann.github.io/DeepDrill/docs/Theory/SmoothColoring.html

The repo also has a Deadline Cloud job template that was used to render the demo video with AWS Deadline Cloud.
## Demo
You can view a demo of the program here:

[![Gameplay](https://img.youtube.com/vi/sHAK_puHJjg/0.jpg)](https://www.youtube.com/watch?v=sHAK_puHJjg)
