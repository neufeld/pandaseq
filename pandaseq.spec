Name:           pandaseq
Version:        2.0
Release:        1%{?dist}
Summary:        Pair-end read assembly
Source:         pandaseq-2.0.tar.bz2
Group:          Applications/Engineering

License:        GPLv3+
URL:            http://github.com/neufeld/pandaseq

BuildRequires:  zlib-devel
BuildRequires:  bzip2-devel
BuildRequires:  libtool-ltdl-devel
BuildRequires:  autoconf
BuildRequires:  automake
BuildRequires:  libtool

%description
PANDA assembles forward and reverse reads from Illumina FASTQ data

%package devel
Summary:        Pair-end read assembly -- Development tools
Requires:       libtool
Requires:	pandaseq

%description devel
PANDA assembles forward and reverse reads from Illumina FASTQ data
This package contains development tools for creating PANDAseq
validation modules. If you are only assembling sequences, this is
not necessary.

%prep
%setup -q
autoreconf -i

%build
%configure
make %{?_smp_mflags}

%install
rm -rf $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT

%files devel
%{_bindir}/pandaxs
%{_includedir}/pandaseq.h
%doc %{_mandir}/man1/pandaxs.1.gz
%doc %{_defaultdocdir}/pandaseq/sample.c
%files
%{_bindir}/pandaseq
%{_libdir}/pandaseq
%doc %{_mandir}/man1/pandaseq.1.gz
%doc %{_defaultdocdir}/pandaseq/README

%changelog
 * Thu Apr 14 2011 Andre Masella <andre@masella.name> 2.0-1
 - First version in C

