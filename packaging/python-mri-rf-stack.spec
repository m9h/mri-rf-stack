Name:           python-mri-rf-stack
Version:        0.1.0
Release:        1%{?dist}
Summary:        MRI RF field analysis and inverse-imaging research tools
License:        MIT
URL:            https://github.com/m9h/mri-rf-stack
Source0:        %{url}/archive/refs/tags/v%{version}.tar.gz#/mri-rf-stack-%{version}.tar.gz
BuildArch:      noarch

BuildRequires:  python3-devel
BuildRequires:  pyproject-rpm-macros
BuildRequires:  python3-pytest

%global _description %{expand:
Research tools for complex MRI RF field analysis, accepted-power normalization,
B1 decomposition, local SAR, JEMRIS coil-map export, and noise-whitened
minimum-norm inverse imaging. Hertzon is an optional separately installed
solver and is not required for field analysis or inverse reconstruction.}

%description %_description

%package -n python3-mri-rf-stack
Summary:        %{summary}
%description -n python3-mri-rf-stack %_description

%prep
%autosetup -n mri-rf-stack-%{version}

%generate_buildrequires
%pyproject_buildrequires

%build
%pyproject_wheel

%install
%pyproject_install
%pyproject_save_files -l mri_rf

%check
%pytest

%files -n python3-mri-rf-stack -f %{pyproject_files}
%doc README.md docs/neurofedora.md docs/validation.md
%{_bindir}/mri-rf

%changelog
* Sun Oct 25 2026 Morgan Hough <morgan.hough@gmail.com> - 0.1.0-1
- Initial analysis and inverse-imaging package; optional solver kept separate
