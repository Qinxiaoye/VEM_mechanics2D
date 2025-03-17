% Corporate article-template for KU Leuven en Kulak
% Stijn Rebry - december 2012
% Andries Vansweevelt - 2018 - aanpassing aan nieuwe huisstijl KU Leuven
% Stijn Rebry - oktober 2021 - afgestoft
% Tips en hints voor een beter resultaat altijd welkom
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{kulakarticle}[2021/10/05 Article class for KU Leuven]
\DeclareOption{kulak}{
   \providecommand{\Logo}{logo}
   }
\DeclareOption{kul}{
   \renewcommand{\Logo}{Humboldt}
   }
% \DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ExecuteOptions{kulak} %Standaardopties
\ProcessOptions\relax
\LoadClass[12pt]{article}
\RequirePackage[a4paper]{geometry}

%\RequirePackage{arev}
%\renewcommand*\familydefault{\sfdefault} %% Only if the base font of the document is to be sans serif
\RequirePackage[T1]{fontenc}
\RequirePackage{hyperref}
\hypersetup{
colorlinks=true,
linkcolor=blue,
citecolor=blue
}
\RequirePackage{titling}
\RequirePackage{multicol}
\RequirePackage{graphicx,xcolor}
	\definecolor{kulgroen}{rgb}{0.1137,0.5529,0.6902} %Oude kleur
	\definecolor{kullichtblauw}{RGB}{106,194,238}
	\definecolor{grijs}{rgb}{0.25,0.25,0.25}
	\definecolor{KULblauw1}{RGB}{29,141,176} %Hoofdblauw nieuwe huisstijl
	\definecolor{KULblauw2}{RGB}{220,231,240} %Secundair blauw: heel lichtblauw
	\definecolor{KULblauw3a}{RGB}{82,189,236} %Tertiair blauw (accent): licht
	\definecolor{KULblauw3b}{RGB}{30,110,135} %Tertiair blauw (accent): donker
	\definecolor{KULblauw3c}{RGB}{12,73,92} %Tertiair blauw (accent): donkerst
\RequirePackage{lastpage}
\RequirePackage{bm}  
\RequirePackage{mathrsfs}
\RequirePackage{amsfonts,amssymb}
\RequirePackage{tikz}
	\usetikzlibrary{arrows} 
	\usetikzlibrary{decorations.pathmorphing}
	\usetikzlibrary{backgrounds}
	\usetikzlibrary{fit}
	\usetikzlibrary{shadows}
	\usetikzlibrary{positioning}
\makeatletter
%\newcommand\faculty[1]{\renewcommand\@faculty{#1}}
%\newcommand\@faculty{}
%\newcommand{\thefaculty}{\@faculty{}}
%\newcommand\group[1]{\renewcommand\@group{#1}}
%\newcommand\@group{}
%\newcommand{\thegroup}{\@group{}}
%\newcommand\subtitle[1]{\renewcommand\@subtitle{#1}}
%\newcommand\@subtitle{}
%\newcommand{\thesubtitle}{\@subtitle{}}
%\newcommand\institute[1]{\renewcommand\@institute{#1}}
%\newcommand\@institute{}
%\newcommand{\theinstitute}{\@institute{}}
\newcommand\address[1]{\renewcommand\@address{#1}}
\newcommand\@address{}
\newcommand{\theaddress}{\@address{}}
\newcommand\location[1]{\renewcommand\@location{#1}}
\newcommand\@location{}
\newcommand{\thelocation}{\@location{}}
%\newcommand\emailaddress[1]{\renewcommand\@emailaddress{#1}}
%\newcommand\@emailaddress{}
%\newcommand{\theemailaddress}{\@emailaddress{}}
\makeatother
%\RequirePackage[lmargin=3.0cm, rmargin=2.cm,tmargin=3cm,bmargin=1.50cm]{geometry}
\setlength{\topmargin}{-1in}
\addtolength{\topmargin}{1.5cm}
\setlength{\oddsidemargin}{-1in}
\addtolength{\oddsidemargin}{3cm}
\setlength{\textwidth}{16cm}
\setlength{\textheight}{23.5cm}
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}
\setlength{\headheight}{1.5cm}
\lhead{}
\rhead{}
\chead{}
\cfoot{}
\rfoot{\footnotesize{\thepage/\pageref{LastPage}}}
% \rhead{\begin{tikzpicture}[inner sep=0mm]
%          \fill[KULblauw3a] (-0.15,0) rectangle (0,1.2);
%          \node[anchor=west] at ( 0.5,0.6) {\includegraphics[height=1.5cm]{Humboldt}};
%       \end{tikzpicture}}
% \lfoot{\begin{tabular}[b]{@{}l@{}}
%     \footnotesize{\thetitle}\\
%     \footnotesize{\theauthor}
%   \end{tabular}}
% \cfoot{}
% \rfoot{\footnotesize{\thepage/\pageref{LastPage}}}

\renewcommand{\maketitle}{
\begin{tikzpicture}[remember picture,overlay]
     \node at (current page.south west)
     {
     \begin{tikzpicture}[remember picture, overlay,inner sep=0mm]
         \fill[KULblauw3a] (17.5,26.87) rectangle (17.65,28.07);
         \node[anchor=north west] at ( 3.0,28.2) {\includegraphics[height=1.5cm]{\Logo}};
         \node[anchor=north west] at ( 8.3,28.2) {\includegraphics[height=1.5cm]{IKM}};
         \node[anchor=north east] at ( 19,28.2) {\includegraphics[height=1.5cm]{Humboldt}};
         \node[anchor=east] at (16.5,27.47) {
         \footnotesize\begin{tabular}{@{}r@{}}\theaddress\end{tabular}
         };
         \node[anchor=south east] at (19,24) {
                   \begin{tabular}{@{}r@{}}
                        \large\bf\thetitle\\[.5cm]
                        \large\theauthor\vphantom{Ky}
                   \end{tabular}};
         \node[anchor=south west] at (3,24) {
					 \begin{tabular}{@{}l@{}}
									\large \thedate\vphantom{Ky}
							   \end{tabular}};
         \draw  (3,23.7) -- ( 19,23.7);
      \end{tikzpicture}
      };
\end{tikzpicture}
\vspace{1cm}\\
{}
% \thispagestyle{empty}
}