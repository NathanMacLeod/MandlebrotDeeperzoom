 COMPILER = g++
 FLAGS = -std=c++17 -lm -pthread -D"NO_INTERACTIVE"
 
 mandlebrot : animation_render.o argument_parsing.o big_float.o main.o mandlebrot.o
 	$(COMPILER) $(FLAGS) -o mandlebrot animation_render.o argument_parsing.o big_float.o main.o mandlebrot.o
 
 .PHONY : clean
 clean : 
 	rm -f *.o mandlebrot
 
 animation_render.o : src/AnimationRender.cpp include/AnimationRender.h
 	$(COMPILER) $(FLAGS) -c -o animation_render.o src/AnimationRender.cpp
 
 argument_parsing.o : src/ArgumentParsing.cpp include/ArgumentParsing.h
 	$(COMPILER) $(FLAGS) -c -o argument_parsing.o src/ArgumentParsing.cpp
 
 big_float.o : src/BigFloat.cpp include/BigFloat.h
 	$(COMPILER) $(FLAGS) -c -o big_float.o src/BigFloat.cpp
 
 main.o : src/main.cpp include/Mandlebrot.h include/InteractiveExplorer.h include/AnimationRender.h include/ArgumentParsing.h
 	$(COMPILER) $(FLAGS) -c -o main.o src/main.cpp
 
 mandlebrot.o : src/Mandlebrot.cpp include/Mandlebrot.h include/stb_image_write.h
 	$(COMPILER) $(FLAGS) -c -o mandlebrot.o src/Mandlebrot.cpp