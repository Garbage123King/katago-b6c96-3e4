```
gcc -g main.c checkresult.c -o main -lz -Wl,--stack,67108864

./main
```



If you want to run 3D demo:
1. Put raylib.h and raylib.dll here;
2. Run command:
```
gcc demo.c -o demo.exe -L. -lraylib -lopengl32 -lgdi32 -lwinmm -mwindows
./demo
```
