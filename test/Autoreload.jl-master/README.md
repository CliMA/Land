Autoreload
===============

Autoreload is a package for autoreloading modules in IJulia. It is intended to allow a workflow where you develop Julia source in some editor, while interacting with it via an IJulia notebook or the command-line REPL. It helps get around annoying type redefinition issues. It is modeled after IPython's autoreload extension.

Installation
=============
In a Julia session, type

```julia
Pkg.update()
Pkg.clone("Autoreload")
```

Usage
=======
First load the package:

```julia
using Autoreload
```

You can then use the ```arequire(filename)```  where you normally would have used ```require``` or ```import```. If you then call ```areload()```, all files included with ```arequire``` will be reloaded if the source files have been modified since their last reload. 

You can use ```smart_reload(filename)``` to reload ```filename``` in a way that avoids type redefinition issues. ```smart_reload``` is automatically called by ```areload```, but you can use it even for files and packages you are not auto-reloading.

A list of files marked for autoreloading can be seen by calling ```arequire()```. A file can be deleted from the autoreload list by calling ```arequire(filename, :off)```.

Module dependencies
====================
There is basic support for handling dependencies between files which are to be reloaded. For example, if a file M3.jl should be loaded only after loading files M1.jl and M2.jl (for example, if M3 imports M1 and M2), you can write

```
arequire("M3", depends_on=["M1", "M2"])
```

M3 will then be auto-reloaded if either M1.jl, M2.jl, or M3.jl is edited, will all three files being reloaded in the correct order.  If an autoreloaded file has ```include``` statements, any file it includes will automatically be determined to be a dependency. This makes it convenient to interactively write a package by calling ```arequire``` with the package name and including the rest of the package files with ```include``` statements in the main package source file.


IJulia integration
===================
If you are using IJulia (recommended), then ```areload()``` will automatically be called before each cell is executed. This behavioral can be toggled on and off by calling ```areload(:on)``` and ```areload(:off)```.

Example
========
In a file called M.jl, I have

```julia
x="First version"
```

In an interactive session, I type

```julia
using Autoreload
arequire("M.jl")
x
```
This will evaluate to "First version".

I then edit M.jl to be

```julia
x="Second version"
```

Then in the same interactive session, I type

```julia
areload()
x
```

and get back "Second version". If I had been using IJulia, the call to ```areload()``` would have been unnecessary.

Package handling
==================
Say you are creating a package organized on disk as ~/.julia/MyPackage/src/[source files].jl. One of the source files will be called MyPackage.jl and is typically loaded to load the rest of the package. If Autoreload finds a file called ``src/MyPackage_reload.jl``, however, then when reloading the package, the package will be reloaded via ``MyPackage_reload.jl`` instead of MyPackage.jl. This allows you to only define constants and code in MyPackage.jl, while MyPackage_reload.jl only reloads code. This allows you to specify custom logic for reloading your package. Here is an example

```julia
arequire("MyPackage") # ~/.julia/MyPackage/src/MyPackage.jl is loaded
...
# make an edit to some file included in MyPackage.jl
areload() # ~/.julia/MyPackage/src/MyPackage_reload.jl is executed instead of ~/.julia/MyPackage/src/MyPackage.jl 
```

This behavior can be disabled by running ```aoptions_set(constants=true)```.


Smart reloading to avoid type redefinition errors
=============================================
Autoreload provides a function called ``smart_reload``. It has similar semantics to ``reload``, but avoids some common issues that make ``reload`` inconvenient for interactive development. 

If you try to reload a type that is already defined in the global scope (e.g, you are auto-reloading a file that defines a type not wrapped in a module), you would normally get an error about redefining a constant. Autoreload will automatically remove the type declaration before reloading your script it is identical to a type that is already defined, avoiding an error.


If you reload a module that defines types, then those type definitions will be stripped out of the module, and the remaining expressions in the reloaded module will be executed in the context of the old module. That way, variables in the global namespace that has the type of a type defined in the module won't have to be redefined when you reload the module. Here is a clarifying example:

A file called M.jl contains:

```julia
module M
type MyType
  x::Int
end

f(var::MyType) = print("First version")
end
```

Then in an interactive session, I have:

```
using Autoreload
arequire("M")
my_var = M.MyType(5)
M.f(my_var)
```

this will print "First version". Now I edit M.jl, and replace it with

```julia
module M
type MyType
  x::Int
end

f(var::MyType) = print("Second version")
end
```

Then in the interactive section, I write

```julia
areload()
M.f(my_var)
```

This will print "Second version". If you had used ```Base.reload("M.jl")``` instead of reloading via ``smart_reload``, "First version" would have been printed in first case, but second case would have resulted in an error. If a file is marked for autoreloading (via ```arequire``), then whenever that file or any file it includes changes, it will be reloaded with ``smart_reload``. 

Limitations
============
Autoreload.jl uses Julia's built-in ```reload``` command, and as such is subject to various limitations inherent in the current Julia architecture. 
