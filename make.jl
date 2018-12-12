using Documenter

makedocs(
    sitename="ECHO RESONANCE Paper 1",
    pages=[
        "Home" => "index.md",
        "Analysis" => [
            "About this section" => "index.md",
            "Data description" => "notebooks/01-data.md",
            ],
        ]
    )
