module ProjectInfo

open Fake.Core

let project = "BioFSharp.Mz"

let testProject = "tests/BioFSharp.Mz.Tests/BioFSharp.Mz.Tests.fsproj"

let summary = "BioFSharp.Mz - modular computational proteomics."

let solutionFile  = "BioFSharp.Mz.sln"

let configuration = "Release"

// Git configuration (used for publishing documentation in gh-pages branch)
// The profile where the project is posted
let gitOwner = "CSBiology"

let gitName = "BioFSharp.Mz"

let gitHome = sprintf "%s/%s" "https://github.com" gitOwner

let projectRepo = sprintf "%s/%s/%s" "https://github.com" gitOwner gitName

let website = "/BioFSharp.Mz"

let pkgDir = "pkg"

let release = ReleaseNotes.load "RELEASE_NOTES.md"

let stableVersion = SemVer.parse release.NugetVersion

let stableVersionTag = (sprintf "%i.%i.%i" stableVersion.Major stableVersion.Minor stableVersion.Patch )

let mutable prereleaseSuffix = ""

let mutable prereleaseTag = ""

let mutable isPrerelease = false
