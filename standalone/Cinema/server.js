'use strict';
// /run/media/jones/Data/Asteroid

const fs = require('fs');
const sqlite3 = require('sqlite3').verbose();
const express = require('express');
const app = express();

const db = new sqlite3.Database(':memory:');
const databasesPath = process.argv[2];
const port = 8888;

app.get('*', function (req, res) {
    const sql = decodeURIComponent(req.url).slice(1);
    console.log(sql);

    db.all(sql, [], (err, rows) => {
        if (err) {
            res.send('[]');
            return;
        }

        console.log('sending ',rows.length,' rows');
        res.send(JSON.stringify(rows));
    });
});

db.serialize(function() {
    db.run("CREATE TABLE VolumeData (Name TEXT, Ariburst TEXT, AsteroidDiameter TEXT, EntryAngle TEXT, Resolution TEXT, CycleTime INT, CDB_Filepath TEXT)");
    db.run("CREATE TABLE Tsunamies (Name TEXT, Ariburst TEXT, AsteroidDiameter TEXT, EntryAngle TEXT, Resolution TEXT, CycleTime INT, CDB_Filepath TEXT)");

    const productTypes = fs.readdirSync(databasesPath);
    for(let productType of productTypes){
        const csvPath = databasesPath + '/' + productType + '/' + productType + '.csv';
        if(fs.existsSync( csvPath )){
            const csv = fs.readFileSync( csvPath, 'UTF-8' );
            const rows = csv.split('\n');
            let insertString = 'INSERT INTO ' + productType + ' VALUES ';
            for(let row of rows){
                insertString += '("'+row.split(',').join('","')+'"),';
            }
            insertString = insertString.slice(0,insertString.length-1) + ';';
            db.run(insertString);
        }
    }

    app.listen(port, () => console.log('server running on port', port));
});



// let sql = `SELECT * FROM VolumeData WHERE time<10000 ORDER BY time`;
