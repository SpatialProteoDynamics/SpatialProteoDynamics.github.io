const fs = require("fs")

function listCSVs() {
   const list = fs.readdirSync(".")
   return list.filter((file) => file.endsWith(".csv"))
}

function readAllFractions(csvFiles) {
    return csvFiles.map(readFraction)
}

function csvToJson(json, variableName, fileName) {
    fs.writeFileSync(fileName, `var ${variableName} = ${JSON.stringify(json)}`)
}

const csvFiles = listCSVs()
const fractions = readAllFractions(csvFiles)
for (let i = 0; i < fractions.length; i++) {
    const fractionName = csvFiles[i].replace(".csv", "")
    csvToJson(fractions[i], `mydata_${fractionName}`, `./output/mydata_${fractionName}.js`)
}

function readFraction(inputFile) {

    const data = fs.readFileSync(inputFile, 'utf8')
    
    const json = {
        datasets: []
    }

    const rows = data.split("\r\n")
    for (const row of rows) {
        const [genenames, site, multiplicity, ...data] = row.split(",")
        if (!json.labels) {
            json.labels = data
        } else {
            json.datasets.push({
                genenames,
                site,
                multiplicity,
                data
            })
        }

    }
    //console.log(`var ${output} = ${JSON.stringify(json)}`)
    return json
}