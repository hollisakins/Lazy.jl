using Documenter
using DocumenterVitepress

makedocs(;
    sitename = "Lazy.jl",
    authors = "Hollis Akins",
    warnonly = true,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/hollisakins/Lazy.jl",
        devbranch = "main",
        devurl = "dev",
    ),
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Getting Started" => "getting-started.md",
        "User Guide" => [
            "Configuration" => "configuration.md",
            "CLI Reference" => "cli.md",
            "Templates" => "templates.md",
            "Filters" => "filters.md",
            "Output Formats" => "output.md",
        ],
        "Advanced Usage" => "advanced.md",
        "Development" => "development.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/hollisakins/Lazy.jl",
    push_preview = true,
)
