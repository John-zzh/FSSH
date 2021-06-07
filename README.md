# Toy surface hopping project in C++

Project Scope:
- Solve the basic fewest switches surface hopping (FSSH) algorithm (Tully, 1990)
- Reproduce trajectory probabilities for Tully's 1D avoided crossing model

We will be using the `eigen` math library for basic linear algebra operations.

## Getting Started: Setting up Your C++ Environment (on a MacOS)

1. Get a text editor with line numbers
  - [Atom](https://atom.io/) is a really nice one. Once you install it, you can add the [clang-format](https://atom.io/packages/clang-format) package which will automatically reformat your C++ code to look pretty every time you save it.

2. Install the `eigen` math library (version 3.3.9) and `cmake` (version 3.19.8)
  - On MacOS: You can quickly install this if you have MacPorts by running the following command:
  ```
  sudo port install eigen3 cmake
  ```

3. Pull the `toyFSSH` project repo from GitLab
  - First, create a fork of the project under your own personal namespace. Check out [this link](https://docs.gitlab.com/ee/user/project/repository/forking_workflow.html) for instructions.
  - Then clone a copy of the forked repository to your personal repository: `git clone <input-from-SSH-key>`

4. Make a directory named `build` on the same level as the `src` directory. From this `build` directory, run `cmake ../src/` from the command line. You should see the following output:
```
-- The C compiler identification is AppleClang 12.0.5.12050022
-- The CXX compiler identification is AppleClang 12.0.5.12050022
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Configuring done
-- Generating done
-- Build files have been written to: /Users/erm/gitlab/toyfssh_cpp/build
```

5. From the `build` directory, run `make`. You should see the following output:
```
Scanning dependencies of target toyFSSH
[ 50%] Building CXX object CMakeFiles/toyFSSH.dir/main.cpp.o
[100%] Linking CXX executable toyFSSH
[100%] Built target toyFSSH
```

6. You should now have an executable file named `toyFSSH` in the `build` directory. Run `./toyFSSH` from the command line, and you should see the following:
```
This is a toy FSSH program in C++

We are using the Eigen3 linear algebra math library.
Here is a happy little matrix :)
  3  -1
2.5 1.5
```

7. Now you can make any edits you want to the `src/main.cpp` file as you write your code program. After you make changes, save the `src/main.cpp` file, then run `make` in the `build` directory. This will compile and link your newly written code and the changes will be reflected in the `toyFSSH` executable file.
