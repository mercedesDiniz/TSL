# Teoria de Sistemas Lineares (TSL)

Este repositório documenta as atividades, materiais de apoio e anotações referentes à disciplina de **Teoria de Sistemas Lineares**, ministrada no primeiro semestre de 2025.

---

## 📚 Ementa da Disciplina

- Introdução aos sistemas dinâmicos e sistemas de controle
- Revisão de álgebra linear
- Descrição de sistemas dinâmicos contínuos e discretos (função de transferência, variáveis de estado, SISO e MIMO)
- Realização de matrizes função de transferência
- Relação entre pólos e autovalores
- Linearidade e linearização
- Conceito de zeros (caso MIMO)
- Estabilidade: entrada-saída, interna, e equação de Lyapunov
- Controlabilidade e observabilidade
- Representações canônicas
- Realimentação de estados (SISO e MIMO)
- Problemas de regulação, seguimento de referência e rejeição de perturbações
- Margens de estabilidade e desempenho robusto
- Controle LQR
- Observador de estados (ordem completa e reduzida)
- Princípio da separação
- Filtro de Kalman e controle LQG

---

## 🎥 Aulas e Materiais de Apoio

- **Playlist oficial das aulas**:  
  [TSL by Antonio Silveira (YouTube)](https://youtube.com/playlist?list=PL6f7H4c5Rq5hSvPXIPeRxZ_pDH8UacinN&si=GWosW05z9zrrPyVH)

- **Repositório Compartilhado de Materiais**:  
  [Dropbox TSL by Antonio Silveira](https://www.dropbox.com/scl/fo/jfhxgk939bcj9j1bv7h46/AB5XhmTI2bTapTI-T_rwnJY?rlkey=cxq7xpizvo950ml8dhoqop8m8&st=ster25ll&dl=0)

- **Softwares de Simulação**:
  - [MATLAB/Simulink](https://www.mathworks.com/)
  - [SCILAB/xCos](http://www.scilab.org/)

- **Revisões e Complementos**:
  - [Sistemas de Controle I - Prof. Andre Dias](https://youtube.com/playlist?list=PL6pEBit8ov96S98YBKH2sGODs6oaS6udX&si=y_bM_2qsZRVasMTp)
  - [Sistemas de Controle II - Prof. Andre Dias](https://youtube.com/playlist?list=PL6pEBit8ov94jES4bFdqTM0aOZxNgNOql&si=tBstPz-kUhHh3b6I)
  - [Controle Digital - Prof. Andre Dias](https://youtube.com/playlist?list=PL6pEBit8ov94P5irBvwEmQfoAAHc9mh_J&si=yt2ioUkwCy3u-rNH)

---

## 📝 Avaliação

- **Trabalhos e Atividades**: 70%  
  (Relatórios e artigos técnicos)
- **Prova Final**: 30%

---

## 📅 Registro das Aulas

### Aula 01 — *26/03/2025*
- Apresentação da disciplina: ementa, cronograma, bibliografia e sistema de avaliação.

---

### Aula 02 — *31/03/2025*
- Projeto exemplo de sistemas de controle:
  - Modelagem
  - Projeto de controladores
  - Análise no domínio do tempo e da frequência

---

### Aula 03 — *02/04/2025*
- Continuação do estudo de caso: modelagem, projeto e análise de controladores.

---

### Aula 04 — *07/04/2025*
- Continuação do estudo de caso:
  - Modelagem experimental (Identificação de Sistemas)
    - Estimador paramétrico por Mínimos Quadrados
  - Projeto de controlador PID digital a partir do PID contínuo por aproximação

---

### Aula 05 — *09/04/2025*
- Continuação da modelagem e controle:
  - Implementação de controlador PID digital
  - Sintonia por tentativa e erro

---

### Aula 06 — *14/04/2025*
- Continuação da análise do projeto:
  - Identificação de sistemas via mínimos quadrados
  - Implementação e sintonia do PID digital
  - Estabilidade relativa:
    - Resposta em frequência: métodos em malha aberta e fechada
  - Sintonias baseadas em modelo:
    - Alocação de polos
    - Métodos ótimos

---

## 📂 Estrutura do Repositório 

```bash
├── aulas/
│   ├── aula02a06/
|   |   ├── imgs/
|   |   └── notas.md
│   └── ...
├── bibliografia/
│   ├── artigos/
│   └── ...
├── projetos/
│   └── ...
├── provas_e_exercicios
│   ├── 2017-1
|   ├── ...
│   └── 2025-1
├── tools
│   ├── daqduino_v2.2.zip
│   └── ...
├── README.md
