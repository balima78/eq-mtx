# usar o {fusen} para criar um package de apoio à analise
# criar um brangh novo para o package a criar
# executar o comando (NOTA - o projecto apenas pode ter letras e pontos, no nome):
fusen::create_fusen(path = ".",
                    template = "minimal",
                    overwrite = TRUE)

# alterar os ficheiros da pasta dev\
# executa o primeiro chunck do ficheiro '0-dev_history.Rmd' para criar o ficheiro DESCRIPTION
# depois de definidas as novas funções no ficheiro 'flat_minimal.Rmd', executar o ultimo chunck para 'inflar o pacote'
# agora já é possível acrescentar dados ao pacote com a função usethis::use_data()
