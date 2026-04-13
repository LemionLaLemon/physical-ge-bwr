var gateway = `ws://${window.location.hostname}/ws`;
var websocket;

window.addEventListener('load', onload);

function onload(event){
    initWs();
}

function initWs(){
    websocket = new WebSocket(gateway);
    websocket.onmessage = onMessage;

    websocket.onclose = () => {
        setTimeout(initWs, 2000);
    }
}

function onMessage(event){
    let data = JSON.parse(event.data);
    //update ui garbage javascript ughhhh
    document.getElementById("CoreFlow").textContent = (((data.RRCA * 34340) + (data.RRCB * 34340))/1000).toFixed(1);
    document.getElementById("APRM").textContent = (data.APRM).toFixed(2);
   
    Object.entries(data.rods).forEach(([rodName, rodData]) => {
        if (rodData.p % 2 == 0 && rodData.p != 48) {
            document.getElementById(rodName).textContent = String(rodData.p).padStart(2, '0');
        } 
        else if (rodData.p == 48) {
            document.getElementById(rodName).textContent = "**";
        }
        else {
            document.getElementById(rodName).textContent = "--";
        }
    });
}
