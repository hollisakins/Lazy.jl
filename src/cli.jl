"""
Command Line Interface for Lazy.jl

Contains CLI parsing, help system, and user interaction functions.
"""

using TOML

function print_help(io, cmd::String = "main")
    printstyled(io, "Lazy.jl \n", bold = true)
    if cmd == "main"
        println("usage: lazy <command> <options>")
        println("")
        println("  fit <options>       Fit the data")
        println("  list-templates      List the available templates")
        println("  list-filters        List the available filters")
        println("  cache-clear         Clear cached template grids")
        println("  -v, --version       Show version information")
        println("  -h, --help          Show this help message")
    elseif cmd == "fit"
        println("")
        println("  -p, --param     Path to the parameter file")
        println("  -t, -- threads  Number of threads to use")
        println("  -y, --yes       Non-interactive mode (auto-resume, no prompts)")
    elseif cmd == "list-templates"
        println("usage: lazy list-templates")
    end

end

function print_ascii(io)
    nthreads = Threads.nthreads()
    println("     __                        _ __ ")
    println("    / /  ____ _____ __  __    (_) / ")
    println("   / /  / __ `/_  // / / /   / / /  ")
    println("  / /__/ /_/ / / // /_/ /   / / /   ")
    println(" /_____|__,_/ /___|__, (_)_/ /_/    ")
    println("                 /____/ /___/       ")
    if nthreads > 1
        println("v$version ($nthreads threads)")
    else
        println("v$version (1 thread)")
    end
    println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
end

function main(argv::Vector{String})
    
    io = stdout    
    if argv == []
        print_ascii(io)
        check_version()
        print_help(io)
        return 0
    end

    # Parse command line arguments
    while length(argv) > 0
        x = popfirst!(argv)
        
        # Return the help message and exit
        if x == "-h" || x == "--help"
            print_ascii(io)
            check_version()        
            print_help(io)
            return 0
        
        # Return version information
        elseif x == "-v" || x == "--version"
            println("Lazy.jl version $version")
            return 0
        
        # 
        elseif x == "fit"
            print_ascii(io)
            check_version()

            if argv == []
                print_help(io, "fit")
                return 0
            end

            param = nothing
            auto_yes = false

            while length(argv) > 0
                y = popfirst!(argv)
                if y == "-p" || y == "--param"
                    if length(argv) < 1
                        throw(LazyError("expected parameter file argument after `-p`"))
                    end
                    param = popfirst!(argv)
                elseif y == "-y" || y == "--yes"
                    auto_yes = true
                elseif y == "-h" || y == "--help"
                    print_help(io, "fit")
                    return 0
                else
                    throw(LazyError("unrecognized argument: $y"))
                end
            end

            if param === nothing
                throw(LazyError("parameter file not specified. Use -p <param_file>"))
            end
            return fit(param; auto_yes=auto_yes)


        elseif x == "list-templates"
            # List the available templates
            print_ascii(io)
            check_version()
            return list_templates()
            
        elseif x == "list-filters"
            # List the available filters
            print_ascii(io)
            check_version()
            return list_filters()
            
        elseif x == "params"
            return params()
            
        elseif x == "cache-clear"
            # Clear template cache
            print_ascii(io)
            check_version()
            removed_count, freed_gb = clear_cache()
            if removed_count > 0
                println("Cleared $removed_count cached template grids ($(round(freed_gb, digits=2)) GB freed)")
            else
                println("No cached template grids found")
            end
            return 0
            
        else
            # Argument not recognized
            throw(LazyError("unrecognized argument: $x"))
        end
    end
    if param == nothing
        throw(LazyError("parameter file not specified. Use -p <param_file>"))
    end

end

function list_filters()
    """
    List the available filters in the filter_files directory.
    """
    # filter_files = glob("*.dat", filterpath)
    filter_directory = load_filters()
    filter_names = sort(collect(keys(filter_directory)))
    
    println("Available filters:")
    for filter_name in filter_names
        filter_description = filter_directory[filter_name]["description"]
        print(rpad(filter_name, 25))
        println(filter_description)
        if !isfile(filterpath * filter_name)
            throw(LazyError("Filter file $filter_name not found"))
        end
    end
end

function list_templates()
    """
    List the available templates in the template_set directory.
    """
    println("Available template sets:")
    template_directory = load_templates()
    template_set_names = sort(collect(keys(template_directory)))
    for template_set_name in template_set_names
        template_files = template_directory[template_set_name]["files"]
        for template_file in template_files
            if !isfile(templatepath * template_file)
                throw(LazyError("Template file $template_file not found"))
            end
        end
        nfiles = length(template_files)
        println("* $template_set_name ($nfiles)")
    end
end

function params()
    """ 
    Read in the example parameter file and print the parameters.
    """
    param_file = @__DIR__() * "/example_params.toml"
    if !isfile(param_file)
        throw(LazyError("Parameter file $param_file not found"))
    end
    f = readlines(param_file)
    for line in f
        println(line)
    end


end

Base.@main