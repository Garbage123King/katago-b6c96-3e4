```
gcc -g main.c checkresult.c -o main -lz -Wl,--stack,67108864

./main
```



If you want to run 3D demo, make sure `glslc` in your path and run:
```
mkdir build
cd build
cmake -S .. -B . -G "Visual Studio 18 2026" -A x64
```

