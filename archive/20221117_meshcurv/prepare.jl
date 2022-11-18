ENV["JULIA_PKG_USE_CLI_GIT"] = "true"

using Pkg

function prepare_ver(rev::AbstractString)
    Pkg.activate(@__DIR__)
    Pkg.add(url="git@github.com:medyan-dev/MEDYAN.jl.git", rev=rev)
    Pkg.instantiate()
end

function prepare()
    @assert length(ARGS) == 1
    if ARGS[1] == "old"
        prepare_ver("72fa5f1cd5b64d5ea2ab494e1394c3bd6053b567")
    elseif ARGS[1] == "new"
        prepare_ver("52da1a000abda1b97eadc35d7b75403f700696e7")
    else
        error("Unrecognized: $(ARGS[1])")
    end
end

prepare()
